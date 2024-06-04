#pragma once

#include <memory>
#include <vector>

// Forward declarations
class Vert;
class Face;
class Storage;

class Edge : public std::enable_shared_from_this<Edge> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vert> vert1, std::weak_ptr<Vert> vert2);
    void delete_();

public:
    int get_id() const;
    std::weak_ptr<Vert> get_vert1() const;
    std::weak_ptr<Vert> get_vert2() const;
    
    void cascade_delete_from_vert(std::weak_ptr<Vert> vert);
    void connect_face(std::weak_ptr<Face> face);
    void disconnect_face(std::weak_ptr<Face> face);

private:
    void check_self_destruction();

    int id_;
    std::weak_ptr<Vert> vert1_;
    std::weak_ptr<Vert> vert2_;
    std::vector<std::weak_ptr<Face>> faces_;
    std::weak_ptr<Storage> storage_;
};